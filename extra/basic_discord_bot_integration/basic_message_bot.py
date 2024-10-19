# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
Basic discord bot setup.

Use https://discord.com/oauth2/authorize?client_id=installation-token-id&permissions=268503040 (permissions for RWX for basic response)

"""

import discord
from discord.ext import commands
import asyncio

TOKEN = "youra-tokena-go-here"
intents = discord.Intents.default()
intents.message_content = True
intents.members = True
bot = commands.Bot(command_prefix='!', intents=intents)

@bot.event
async def on_ready():
    print(f"Logged in as {bot.user} ({bot.user.id})")
    print("Bot is ready and running!")

@bot.event
async def on_message(message):
    if message.author == bot.user:
        return
    if message.channel.name == 'general':
        print(f"Received message in #general from {message.author.display_name}: {message.content}")
        await message.reply(f"Hello, {message.author.display_name}! You said: {message.content}")
    await bot.process_commands(message)

if __name__ == "__main__":
    async def run_bot():
        try:
            await bot.start(TOKEN)
        except discord.PrivilegedIntentsRequired:
            print("Privileged intents are required and have not been enabled.")
            print("Enable intents at https://discord.com/developers/applications/")
            sys.exit(1)
        except discord.LoginFailure:
            print("Login failed: Invalid token provided.")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)

    try:
        asyncio.run(run_bot())
    except KeyboardInterrupt:
        print("Bot is shutting down...")

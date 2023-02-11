import fmt
import k8s.io.apimachinery.pkg.runtime.schema as schema


class KindVisitor:
    def visit_daemon_set(self, kind):
        pass
    
    def visit_deployment(self, kind):
        pass
    
    def visit_job(self, kind):
        pass
    
    def visit_pod(self, kind):
        pass
    
    def visit_replica_set(self, kind):
        pass
    
    def visit_replication_controller(self, kind):
        pass
    
    def visit_stateful_set(self, kind):
        pass
    
    def visit_cron_job(self, kind):
        pass
